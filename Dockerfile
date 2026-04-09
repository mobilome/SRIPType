FROM python:3.11-slim

LABEL maintainer="SRIPType Team"
LABEL description="SRIPType - A bioinformatics toolkit for sequence typing and analysis"

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        procps \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt/sriptype

# Copy project files
COPY . .

# Install SRIPType
RUN pip install --no-cache-dir . && \
    chmod +x sriptype

# Add to PATH
ENV PATH="/opt/sriptype:${PATH}"

# Default working directory for data
WORKDIR /data

ENTRYPOINT ["sriptype"]
CMD ["--help"]
