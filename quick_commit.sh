if [ $# -eq 0 ]; then  # No arguments
  git commit -am "Quick commit"
else
  git commit -am "$1"
fi
