FROM python:3.8

EXPOSE 443
EXPOSE 5000
EXPOSE 8080

# set base directory
WORKDIR /ADZE

# Environmental variables
ENV FEMM /root/.wine/drive_c/femm42/bin/femm.exe
ENV DISPLAY :127
ENV WINEPREFIX /root/.wine
ENV WINEARCH win32
# this needed when python packages get compiled to save space
ENV CLFAGS '-Os -g0 -Wl, --strip-all -I/usr/include:/usr/local/include -L/usr/lib:/usr/local/lib'

# Basic packages
RUN apt-get update \
    && apt-get install -y \
        wget \
        dpkg \
        xvfb \
        cabextract

# PYTHON
COPY requirements.txt .
# Use this to fast install packages
RUN pip3 install --no-cache-dir -r requirements.txt
# Use this in production. It will compile every packages from scratch. Note: Be
# aware this step can be time consuming.
# RUN pip3 install --no-cache-dir --compile --global-option=build_ext -r requirements.txt
RUN echo "/ADZE/" > /usr/local/lib/python3.8/site-packages/adze.pth


# WINE
RUN dpkg --add-architecture i386 \ 
    && apt-get update \
    && apt-get install -y --no-install-recommends wine wine32 \
    && wget  https://raw.githubusercontent.com/Winetricks/winetricks/master/src/winetricks \
    && chmod +x winetricks \
    && apt install -y cabextract \
    && sh winetricks mfc90 \
    && apt remove -y cabextract && rm -f winetricks \ 
    && apt autoremove -y \
    && apt clean -y \
    && rm -rf /var/lib/apt/lists/*


# Install FEMM
ADD resources/femm42.tar.xz /root/.wine/drive_c/

# Install ADZE-modeler
COPY adze_modeler/ adze_modeler/

# get the size of the python packages
# RUN du -sh /usr/local/lib/python3.8/site-packages/* | sort -h

# TEST
COPY start.sh . 
COPY applications/PowerTransformer/snapshots/dev/P_dev.lua test.lua
RUN chmod +x start.sh
CMD ["sh", "start.sh"]


