FROM python/3.8.4-alpine3.12

RUN pip install -r requirements.txt
RUN pip install .