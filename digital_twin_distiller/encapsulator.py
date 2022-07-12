import json
import os.path
import subprocess
import time
import traceback
from pathlib import Path
from typing import Optional

import uvicorn
from fastapi import FastAPI, Response, status
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Extra, ValidationError

from digital_twin_distiller.modelpaths import ModelDir


class InputJsonSim(BaseModel):
    """
    Class for validating the input sent to the /process_sim endpoint for SimulationProject.
    """

    simulation: Optional[dict] = {"type": "default"}
    model: Optional[dict] = {}
    tolerances: Optional[dict] = {
        "type": "ff",
        "parameters": {},
        "variables": [],
    }
    misc: Optional[dict] = {"processes": 4, "cleanup": True, "exportname": None}
    # TODO: get this version number from poetry, what do we want from this version
    version: Optional[str] = "2021.12"

    # Setting for keeping the additional keys in the input json intact
    class Config:
        extra = Extra.allow


# Defining the API
app = FastAPI(title="{} API", docs_url="/apidocs", redoc_url=None)

tags_metadata = [
    {
        "name": "process_sim",
        "description": "Run simulations for the API.",
        "externalDocs": {
            "description": "Find out more",
            "url": "http://montana.ai",
        },
    },
    {"name": "ping", "description": "Endpoint for pinging server."},
    {
        "name": "docs",
        "description": "Endpoint for OpenAPI documentation.",
        "externalDocs": {
            "description": "Find out more",
            "url": "http://montana.ai",
        },
    },
    {
        "name": "root",
        "description": "Test page for the API. Endpoint called by the main page of the API test page.",
    },
    {"name": "docs", "description": "Endpoint for the project documentation."},
]


@app.post("/process_sim", include_in_schema=True, tags=["process_sim"])
async def process_sim(item: InputJsonSim):
    """
    Endpoint for performing the project.run() method on data sent for the API in JSON format.
    The endpoint performs automatic input validation via the Item class.
    """
    data = json.loads(item.json())
    app.project._output.clear()
    try:
        app.project._input = data
        app.project.update_input()
        app.project.run()
    except Exception as e:
        app.project._output["exception"] = {
            "type": e.__class__.__name__,
            "message": str(e),
            "traceback": traceback.format_exc(),
        }
    finally:
        return app.project._output


@app.get("/ping", include_in_schema=True, tags=["ping"])
def ping():
    """
    Pings the server to check if it is available.
    """
    result_json = {}
    result_json["call_time"] = time.ctime()
    result_json["msg"] = "The API is working."
    return result_json


class Encapsulator:
    """
    Server for running a custom project as an API.
    """

    def __init__(self, project):
        """
        :param project: SimulationProject instance
        """
        self.app = app
        self.app.project = project
        self.app.title = self.app.title.format(project.app_name)
        # self.workers = self.number_of_workers()
        self.host = "127.0.0.1"
        self.port = 5000
        self.cert_file_path = None
        self.key_file_path = None

        self.set_endpoint_for_docs(ModelDir.DOCS)

    def set_cert_file_path(self, cert_file_path):
        self.cert_file_path = cert_file_path

    def set_key_file_path(self, key_file_path):
        self.key_file_path = key_file_path

    def set_endpoint_for_docs(self, docs_path, endpoint="/", build=True):
        """
        Builds mkdocs documentation and deploys the documentation at the given endpoint.
        :param docs_path: path to the docs folder containing mkdocs.yml
        :param endpoint: endpoint to publish built docs to
        :param build: whether to build documentation or not
        :return:
        """
        if build:
            self.build_docs(docs_path)
        site_path = Path(docs_path).joinpath("site")
        self.app.mount(endpoint, StaticFiles(directory=site_path, html=True, check_dir=True), name="documentation")

    def build_docs(self, docs_path):
        """
        Build the documentation with mkdocs.
        """
        cwd = Path(os.getcwd())
        os.chdir(docs_path)
        subprocess.run("mkdocs build", shell=True, check=True)
        os.chdir(cwd)
        (docs_path / "site" / "images").resolve().mkdir(exist_ok=True)

    def set_host(self, host: str):
        """
        Set the IP address of the host.
        :param host: e.g. 127.0.0.1
        """
        self.host = str(host)

    def set_port(self, port: int):
        """
        Set the port.
        :param port: int, e.g. 5000
        """
        self.port = int(port)

    def run(self):
        """
        Running the application that is running the specified input project's run method.
        :return: None
        """
        if self.key_file_path and self.cert_file_path:
            uvicorn.run(
                self.app,
                host=self.host,
                port=self.port,
                log_level="info",
                ssl_keyfile=self.key_file_path,
                ssl_certfile=self.cert_file_path,
            )
        else:
            uvicorn.run(self.app, host=self.host, port=self.port, log_level="info")

    def __call__(self, *args, **kwargs):
        self.run()


if __name__ == "__main__":
    uvicorn.run("server:app", host="127.0.0.1", port=5000, log_level="info")
