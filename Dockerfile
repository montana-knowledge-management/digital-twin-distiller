FROM python:3.8

EXPOSE 443
EXPOSE 5000
EXPOSE 8080

RUN apt update && apt install -y xvfb software-properties-common
COPY requirements.txt /root/
COPY adze_modeler/ /root/ADZE/adze_modeler/

# WINE
RUN dpkg --add-architecture i386
RUN apt-get update
RUN apt install -y --install-recommends wine64 wine32
RUN wine --version

RUN wget  https://raw.githubusercontent.com/Winetricks/winetricks/master/src/winetricks
RUN chmod +x winetricks
RUN apt install -y cabextract
RUN sh winetricks mfc90
RUN apt remove -y cabextract && rm -f winetricks



# Install FEMM
ADD resources/femm42.tar.xz /root/.wine/drive_c/


ENV VIRTUAL_ENV=/root/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN echo "/root/ADZE/" > $VIRTUAL_ENV/lib/python3.8/site-packages/adze.pth
RUN pip3 install -r /root/requirements.txt

COPY applications/PowerTransformer/ /root/model/
WORKDIR /root/model

# X virtual frame buffer to run femm headless
COPY applications/PowerTransformer/start.sh .
RUN chmod +x start.sh

CMD ["./start.sh"]
