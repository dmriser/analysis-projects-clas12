FROM rootproject/root-ubuntu16 

# Run as superuser 
USER root 

ARG username=user
RUN userdel builder && useradd --create-home --home-dir /home/${username} ${username}
ENV home /home/${username}

# Switch to new user 
USER ${username}
WORKDIR /data

# CMD ["python", "pres.py"]
CMD ["python", "test.py"]