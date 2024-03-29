FROM jfreeman88/saher_collab:0.1


# copy the current directory contents into the container at /app
COPY . /app

# Set the working directory
WORKDIR /app
