FROM python:3.8

EXPOSE 443
EXPOSE 5000
EXPOSE 8080

# set base directory
WORKDIR /ADZE

# Basic packages
RUN apt-get update \
    && apt-get install -y \
        wget \
        dpkg \
        xvfb \
        cabextract

# WINE
RUN dpkg --add-architecture i386 \ 
    && apt-get update \
    && apt-get install -y --install-recommends wine64 wine32 \
    && wget  https://raw.githubusercontent.com/Winetricks/winetricks/master/src/winetricks \
    && chmod +x winetricks \
    && apt install -y cabextract \
    && sh winetricks mfc90 \
    && apt remove -y cabextract && rm -f winetricks \ 
    && apt autoremove -y \
    && apt clean -y \
    && rm -rf /var/lib/apt/lists/*

# Python
COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt
RUN echo "/ADZE/" > /usr/local/lib/python3.8/site-packages/adze.pth

# Install FEMM
ADD resources/femm42.tar.xz /root/.wine/drive_c/

# Install ADZE-modeler
COPY adze_modeler/ adze_modeler/

ENV FEMM /root/.wine/drive_c/femm42/bin/femm.exe
ENV DISPLAY :127

# TEST
COPY start.sh . 
COPY applications/PowerTransformer/snapshots/dev/P_dev.lua test.lua
RUN chmod +x start.sh
CMD ["sh", "start.sh"]
