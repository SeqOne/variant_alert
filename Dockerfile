FROM python:3.6-stretch
LABEL maintainer='raphael.lanos@seqone.fr'

RUN pip install poetry

WORKDIR /root

COPY . .
RUN poetry install

ENTRYPOINT ["poetry", "run"]
