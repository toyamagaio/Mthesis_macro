#! /bin/csh
sed -i -e "1s/^/#/g" ${1}
sed -i -e "s/,/ /g" ${1}
sed -i -e "s/inf/9999/g" ${1}
