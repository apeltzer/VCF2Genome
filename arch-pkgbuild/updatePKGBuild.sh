#!/bin/bash
#Replace version information appropriately in PKGBUILD and starterScript
VERSION=$(grep -e '^version' ../build.gradle | tr -d "'" | tr -d "version ")
VERSIONREL=1
#PKGBUILD
sed -i "s/^pkgver=.*/pkgver=$VERSION/" PKGBUILD
sed -i "s/^pkgrel=.*/pkgrel=$VERSIONREL/" PKGBUILD
cp ../build/libs/VCF2Genome*.jar . 
makepkg
