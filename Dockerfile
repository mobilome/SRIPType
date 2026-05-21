FROM python:3.11-slim

LABEL maintainer="Mobilome Team"
LABEL description="scgt - SINE-RIP chip genotyper, a bioinformatics toolkit for 51k-SINE-RIP chip genotyping"

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        procps \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt/scgt

# Copy project files
COPY . .

# Install scgt
RUN pip install --no-cache-dir . && \
    chmod +x scgt

# Add to PATH
ENV PATH="/opt/scgt:${PATH}"

# Default working directory for data
WORKDIR /data

ENTRYPOINT ["scgt"]
CMD ["--help"]
