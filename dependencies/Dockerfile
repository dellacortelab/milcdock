FROM nvidia/cuda:11.2.1-runtime-ubuntu20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN \
  # Update package list
  apt-get update -y && \
  # Install...
  apt-get install -y \
  git \
  vim \
  wget \
  python3 \
  python3-dev \
  python3-pip \
  tmux \
  # Remove package lists
  && rm -rf /var/lib/apt/lists/*

# Install conda/mamba
RUN wget -O Mambaforge.sh  "https://github.com/conda-forge/miniforge/releases/download/22.11.1-1/Mambaforge-22.11.1-1-$(uname)-$(uname -m).sh"
RUN bash Mambaforge.sh -b -p $HOME/conda
RUN $HOME/conda/bin/conda init bash

COPY milcdock_env.yml /tmp/milcdock_env.yml
RUN $HOME/conda/bin/mamba env create -f /tmp/milcdock_env.yml

RUN echo "export PATH=\$HOME/conda/bin:\$PATH" >> $HOME/.bashrc
RUN echo "conda init bash" >> $HOME/.bashrc
RUN echo "conda activate milcdock" >> $HOME/.bashrc

# Install vina
RUN wget https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz
RUN tar -xf autodock_vina_1_1_2_linux_x86.tgz
RUN mv /autodock_vina_1_1_2_linux_x86/bin/vina /usr/local/bin/vina

# Install ledock
RUN wget http://www.lephar.com/download/ledock_linux_x86
RUN wget http://www.lephar.com/download/lepro_linux_x86
RUN chmod +x ledock_linux_x86
RUN chmod +x lepro_linux_x86
RUN mv ledock_linux_x86 /usr/local/bin/ledock
RUN mv lepro_linux_x86 /usr/local/bin/lepro

# Install autodock
RUN wget https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/autodocksuite-4.2.6-x86_64Linux2.tar
RUN tar -xf autodocksuite-4.2.6-x86_64Linux2.tar
RUN mv x86_64Linux2/* /usr/local/bin/

# Install PLANTS
RUN wget http://www.tcd.uni-konstanz.de/plants_download/download/PLANTS1.2_64bit
RUN mv PLANTS1.2_64bit /usr/local/bin/plants
RUN chmod +x /usr/local/bin/plants

# Install MGTools (necessary for running vina and autodock)
RUN wget https://ccsb.scripps.edu/mgltools/download/495/MGLTools-1.5.6-Linux-x86_64.tar.gz
RUN tar -xf MGLTools-1.5.6-Linux-x86_64.tar.gz
RUN cd mgltools_x86_64Linux2_1.5.6 && ./install.sh
RUN echo "export PATH=\$PATH:/mgltools_x86_64Linux2_1.5.6/bin" >> $HOME/.bashrc
RUN echo "export PATH=\$PATH:/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24" >> $HOME/.bashrc

# rDock is installed with conda


ENTRYPOINT ["/bin/bash"] 