# VERSION 1.0, tagged with R version 3.5.2.
# AUTHORS: NIH Metabolomics Data Portal Hackathon Team
# DESCRIPTION: Shiny container with packages needed to perform metabaolmics analysis and visualziation
# BUILD: docker build --rm -t gautham1/shiny-metabolomics:3.5.2 .

FROM rocker/shiny:3.5.2

RUN R -e "install.packages(c('data.table','DT','ggplot2','grid','Hmisc','igraph','remotes','shiny','shinydashboard'))"
RUN R -e "remotes::install_github('NCBI-Hackathons/Metabolomics-Data-Portal')"

RUN localedef -i en_US -f UTF-8 en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
RUN mkdir -p /var/lib/shiny-server/bookmarks && \
  chown shiny:0 /var/lib/shiny-server/bookmarks && \
  chmod g+wrX /var/lib/shiny-server/bookmarks && \
  mkdir -p /var/log/shiny-server && \
  chown shiny:0 /var/log/shiny-server && \
  chmod g+wrX /var/log/shiny-server


EXPOSE 3838

COPY shiny-server.sh /usr/bin/shiny-server.sh

USER shiny

CMD ["/usr/bin/shiny-server.sh"]