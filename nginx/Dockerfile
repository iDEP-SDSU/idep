# Nginx with sticky module
# More detail about sticky module: https://bitbucket.org/nginx-goodies/nginx-sticky-module-ng/overview
FROM nginx:1.15-alpine

# File Author / Maintainer
MAINTAINER Ge Lab

# Copy custom configuration file from the current directory
COPY nginx.conf /etc/nginx/nginx.conf

# for HTTPS
RUN mkdir -p /etc/pki/nginx/private
COPY fullchain.pem /etc/pki/nginx/server.crt
COPY privkey.pem /etc/pki/nginx/private/server.key
RUN chown nginx /etc/pki/nginx/server.crt /etc/pki/nginx/private/server.key