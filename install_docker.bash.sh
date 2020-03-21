#!/bin/bash

sudo addgroup --system docker

sudo adduser ${USER} docker

newgrp docker

sudo snap install docker



docker run hello-world