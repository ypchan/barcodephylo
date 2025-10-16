ls scripts | grep '.py' | while read s; do
    rm ~/bin/${s}
    chmod 755 scripts/${s}
    echo "chmod 755 scripts/${s}"
    ln -s $(pwd)/scripts/${s} ~/bin/${s}
    echo "ln -s $(pwd)/scripts/${s} ~/bin/${s}"
done

chmod 755 barcodephylo.sh
echo "chmod 755 barcodephylo.sh"
ln -s $(pwd)/barcodephylo.sh ~/bin/barcodephylo.sh
echo "ln -s $(pwd)/barcodephylo.sh ~/bin/barcodephylo.sh"