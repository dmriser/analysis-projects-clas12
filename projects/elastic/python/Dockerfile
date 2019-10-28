FROM rootproject/root-ubuntu16 

# Run as superuser 
USER root 

ARG username=user
RUN userdel builder && useradd --create-home --home-dir /home/${username} ${username}
ENV home /home/${username}

# Switch to new user 
USER ${username}
WORKDIR /data

CMD ["python", "mon.py", "-i=rga-in.root", "-o=rga-in"]
# CMD ["python", "mon.py", "-i=sim0.root", "-o=sim0"]
# CMD ["python", "test.py"]