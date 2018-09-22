git add .
if [ $# -eq 0 ]  # No arguments
then
  git commit -m "Quick commit"
else
  git commit -m "$1"
fi
