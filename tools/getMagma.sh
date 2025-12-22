#!/usr/bin/bash
linuxURL=https://vu.data.surf.nl/public.php/dav/files/lxDgt2dNdNr6DYt/?accept=zip
macosURL=https://vu.data.surf.nl/public.php/dav/files/TOH4SuvczAKE29d/?accept=zip
windowsURL=https://vu.data.surf.nl/public.php/dav/files/TOH4SuvczAKE29d/?accept=zip

# Download Magma based on OS
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    wget -c -O magma.zip $linuxURL
    ls -lh magma.zip
elif [[ "$OSTYPE" == "darwin"* ]]; then
    wget -O magma.zip $macosURL
elif [[ "$OSTYPE" == "cygwin" || "$OSTYPE" == "msys" || "$OSTYPE" == "win32" ]]; then
    wget -O magma.zip $windowsURL
else
    echo "Unsupported OS type: $OSTYPE"
    exit 1
fi
unzip -f magma.zip -d inst/bin/magma