if [ $# -eq 0 ]  # No arguments
then
  git commit -am "Quick commit"
else
  git commit -am "$1"
fi
