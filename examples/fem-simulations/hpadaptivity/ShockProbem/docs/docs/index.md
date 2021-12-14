# Simple Api for Adze Projects

The goal of this API is to provide a simple and easy to use interface, which
can help to automatically deploy ADZE-models as a  digital twin. The server
waits for a single json, which contains the input, then gives back the
computation results.

![tutorial](images/tutorial.gif)

## Endpoints

* `/` endpoint, which contains the project documentation, by default, this page
  contains the documentation of the server module. This can be easily replaced
  by a static html page. ADZE-models supports MkDocs (https://www.mkdocs.org/)
  documentation, however, other documentation generator based outputs can be
  inserted to this endpoint.

* `/process` endpoint, where the application waits for an input json as a
  single input.

* `/apidocs` endpoint, which describes the roles and the basic functions of the
  integrated api endpoint as an OpenAPI documentation. Moreover, this html site
  contains a small, integrated **test interface**.

* `/ping` endpoint, to test the functionality and the accessibility of the
  realized server from the client side.

### ADZE-Server stands on the shoulders of giants:

    FastAPI for the web parts.
    Pydantic for the data parts and the JsonSchma validation.
    OpenAPI for the API documentation and testing the basic functionality.

#### Documentation

`set_project_mkdocs_dir_path(path)`:  the generated documentation should be
placed in this directory, the function waits an MkDocs Project which has a
following structure:

      mkdocs.yml      # The configuration file.
        site/         # The generated project files
          assets/
          search/
          index.html  # The endpoint shows this page as the documentations homepage
          ...         # Other generated html pages
        docs/
          index.md    # The documentation homepage.
          ...         # Other markdown pages, images and other files.

`build_docs`: Builds the MkDocs documentation automatically.


#### Secure all endpoints with 'https' communication:

* server.set_key_file_path(key.pem)
* server.set_cert_file_path(files(cert.pem)
* server.set_json_validator(NewJson)
