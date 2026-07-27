# OCI index digest resolved from Docker Hub on 2026-07-28.  The tag documents
# the intended CPython release; the digest prevents later tag movement.
FROM python:3.11.9-slim-bookworm@sha256:8fb099199b9f2d70342674bd9dbccd3ed03a258f26bbd1d556822c6dfc60c317

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONHASHSEED=0 \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 \
    BLIS_NUM_THREADS=1

WORKDIR /workspace

COPY requirements/constraints-container.txt /tmp/constraints-container.txt
RUN python -m pip install --no-cache-dir "pip==26.1.2"
ENV PIP_CONSTRAINT=/tmp/constraints-container.txt

# Copy every input covered by source_fingerprint(), plus the tests.
COPY pyproject.toml ./
COPY Dockerfile ./
COPY .dockerignore ./
COPY requirements ./requirements
COPY src ./src
COPY experiments ./experiments
COPY configs ./configs
COPY scripts ./scripts
COPY tests ./tests

RUN python -m pip install --no-cache-dir ".[dev]"

CMD ["python", "-m", "pytest"]
