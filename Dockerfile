FROM python:3.9

EXPOSE 443
EXPOSE 5000

COPY requirements.txt ./
RUN apt-get update &&  apt-get install -y default-jre
RUN apt-get install -y git git-lfs

RUN pip install -r requirements.txt

COPY docs/ docs/
COPY project/ project/
COPY extractors/ extractors/
COPY fixers/ fixers/
COPY classifiers/ classifiers/
COPY preprocessors/ preprocessors/
COPY resources/ resources/
WORKDIR /project

CMD ["python", "wk_labeler_project.py"]
