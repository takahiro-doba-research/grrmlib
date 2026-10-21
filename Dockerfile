FROM python:3.13-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /work

COPY requirements.txt .
RUN pip install -r requirements.txt

# grrmlib 本体はマウントしたソースを PYTHONPATH で読む（編集がすぐ反映される）
ENV PYTHONPATH=/work/src
