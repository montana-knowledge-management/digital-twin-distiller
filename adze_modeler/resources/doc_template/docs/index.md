# Simple Api for Adze Projects

The goal of this API is to provide a simple and easy to use interface, which can help to automatically deploy
Distiller Project class based machine learning models as a production ready application. The server waits for a single
json, which contains the input text under a `Text` key, then the application reads the output from the output_stack of
the Distiller Project.

## Basic Usage

We can run a simple ExampleProject, which based on the `AbstractProject` classes in the following way:

    from distiller.server import Server
    from example_project import ExampleProject

    server = Server(ExampleProject)
    server.run()

## Endpoints

* `/` endpoint, which contains the project documentation, by default, this page contains the documentation of the server
  module. This can be easily replaced by a static html page. Distiller project supports MkDocs (https://www.mkdocs.org/)
  documentation, however, other documentation generator based outputs can be inserted to this endpoint.

* `/process` endpoint, where the application waits for an input json as a single input for the machine learning based text
  processing task and gives back the distilled text. This `/process` endpoint waits for a json, which contains a `Text` key.

* `/apidoc` endpoint, which describes the roles and the basic functions of the integrated api endpoint as an OpenAPI
  documentation. Moreover, this html site contains a small, integrated **test interface**.

* `/ping` endpoint, to test the functionality and the accessibility of the realized server from the client side.

## Usage of the Integrated Test Interface in `\apidoc`



### Distiller Server stands on the shoulders of giants:

    FastAPI for the web parts.
    Pydantic for the data parts and the JsonSchma validation.
    OpenAPI for the API documentation and testing the basic functionality.

## Customization

Using the integrated methods and parameters of the **Server** class, it can be easily customized:

    EXAMPLE I

    TITLE/HOST/PORT, etc...

#### Documentation

`set_project_mkdocs_dir_path(path)`:  the generated documentation should be placed in this directory, the function waits
an MkDocs Project which has a following structure:

      mkdocs.yml      # The configuration file.
        site/         # The generated project files
          assets/
          search/
          index.html  # The endpoint shows this page as the documentations homepage
          ...         # Other generated html pages
        docs/
          index.md    # The documentation homepage.
          ...         # Other markdown pages, images and other files.

#### Secure all endpoints with 'https' communication:

* server.set_key_file_path(key.pem)
* server.set_cert_file_path(files(cert.pem)
* server.set_json_validator(NewJson)

    EXAMPLE II
