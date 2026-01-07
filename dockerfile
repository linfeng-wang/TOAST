# Base Python image
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Copy source code and metadata
COPY toast_amplicon/ toast_amplicon/
COPY pyproject.toml .
COPY PYPI.md .
COPY LICENSE .

# Install flit and dependencies
RUN pip install --upgrade pip \
 && pip install .

# If you have a CLI script installed by flit, this is the default command
CMD ["toast", "--help"]
