# The container created from this image is used as a jump point within the
# kubernetes cluster to submit various ad-hoc jobs.

FROM nextflow/nextflow:23.10.1

# Adding a few command line utilities
RUN yum update -y && yum install -y --setopt=skip_missing_names_on_install=False \
nano-2.9.8-2.amzn2.0.1.x86_64 \
tar-2:1.26-35.amzn2.0.3.x86_64 \
less-458-9.amzn2.0.4.x86_64 \
wget-1.14-18.amzn2.1.x86_64 \
unzip-6.0-57.amzn2.0.1.x86_64 \
zip-3.0-11.amzn2.0.2.x86_64 \
rclone-1.55.1-1.amzn2.0.1.x86_64

# Installing the aws cli. Using version >= 2.13.0 so that environment variable AWS_ENDPOINT_URL is supported
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.13.0.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install