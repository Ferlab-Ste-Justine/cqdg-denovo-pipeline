# The container created from this image is used as a jump point within the
# kubernetes cluster to submit various ad-hoc jobs.

FROM nextflow/nextflow:23.10.1

# Adding a few command line utilities
RUN yum update -y && yum install -y nano-2.9.8-2.amzn2.0.1.x86_64 \
tar-2-1.26-35.amzn2.0.3.x86_64 \
less.x86-458-9.amzn2.0.3.x86_64