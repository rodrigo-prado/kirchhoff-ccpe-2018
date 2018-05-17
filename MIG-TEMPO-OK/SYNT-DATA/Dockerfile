# Use an official Python runtime as a parent image
FROM ubuntu:14.04
#FROM phusion/baseimage:0.9.17

# Set the working directory to /app
WORKDIR /kirchhoff

# Copy the current directory contents into the container at /app
ADD . /kirchhoff

# Install any needed packages specified in requirements.txt
RUN apt-get update && apt-get install -y \
	python-pip curl apt-transport-https time openmpi-bin openmpi-common libopenmpi-dev openssh-client openssh-server python-dev vim \
	#&& pip install --trusted-host pypi.python.org -r requirements.txt \
	&& pip install -r requirements.txt \
	&& curl -o /root/GPG-PUB-KEY-INTEL-PSXE-RUNTIME-2018 -SL "https://apt.repos.intel.com/2018/GPG-PUB-KEY-INTEL-PSXE-RUNTIME-2018" \
	&& apt-key add /root/GPG-PUB-KEY-INTEL-PSXE-RUNTIME-2018 \
	&& /bin/bash -c 'echo deb https://apt.repos.intel.com/2018 intel-psxe-runtime main > /etc/apt/sources.list.d/intel-psxe-runtime-2018.list \
	&& echo deb [arch=all] https://apt.repos.intel.com/2018 intel-psxe-runtime main >> /etc/apt/sources.list.d/intel-psxe-runtime-2018.list' \
	&& apt-get update \
	&& apt-get install -y intel-psxe-runtime \
	&& /bin/bash -c 'source /opt/intel/psxe_runtime/linux/bin/psxevars.sh'

#RUN curl -o /root/GPG-PUB-KEY-INTEL-PSXE-RUNTIME-2018 -SL "https://apt.repos.intel.com/2018/GPG-PUB-KEY-INTEL-PSXE-RUNTIME-2018" \
#	&& apt-key add /root/GPG-PUB-KEY-INTEL-PSXE-RUNTIME-2018 \
#	&& /bin/bash -c 'echo deb https://apt.repos.intel.com/2018 intel-psxe-runtime main > /etc/apt/sources.list.d/intel-psxe-runtime-2018.list \
#	&& echo deb [arch=all] https://apt.repos.intel.com/2018 intel-psxe-runtime main >> /etc/apt/sources.list.d/intel-psxe-runtime-2018.list' \
#	&& apt-get update \
#	&& apt-get install -y intel-psxe-runtime 

# Make port 80 available to the world outside this container
# EXPOSE 80

# Define environment variable
# ENV NAME World

# Run app.py when the container launches
#CMD ["run-openmp-auto.sh", "8", "5", "2"]
# CMD ["python", "run-batch.py"]
