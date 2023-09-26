# AMAA
## Download and install Docker
Docker is available at https://docker.com
## Get AMAA Docker image
    docker image pull viva1109/amaa:1.02
## Get additional tools and example datsets from github
    git clone https://github.com/viva1109/amaa
## Docker run
    docker run --name amaa -v $(pwd):/home/ -it viva1109/amaa
