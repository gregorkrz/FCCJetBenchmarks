FROM cern/alma9-base

RUN dnf -y install --nodocs --setopt=install_weak_deps=False \
        mesa-libGLU \
        make \
        cmake \
        gcc \
        gcc-c++ \
        glibc-devel \
        glibc-headers \
        kernel-headers \
        libstdc++-devel && \
    dnf clean all && \
    rm -rf /var/cache/dnf /var/log/dnf* /var/log/yum.*


# Set default command (optional)
CMD ["/bin/bash"]

