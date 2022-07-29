
FROM python:3.9.13

# set up working directory
WORKDIR /usr/src/app

# install the digital-twin-distiller package
RUN pip install digital-twin-distiller

# copy model files into the working directory
COPY . .

# start the server
CMD [ "python", "./simulation.py" ]
