FROM python:3.12

COPY app.py requirements.txt gunicorn_config.py /app/

RUN pip install -r /app/requirements.txt

WORKDIR /app
CMD gunicorn --worker-tmp-dir /dev/shm --config /app/gunicorn_config.py app:app

