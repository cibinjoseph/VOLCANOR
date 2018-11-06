if [ $# -eq 0 ]; then      # No arguments
  git commit -am "Quick commit"
elif [[ $1 == '-t' ]] | [[ $1 == 't' ]]; then
  git commit -am 'testcase'
elif [[ $1 == '-a' ]] | [[ $1 == 'a' ]]; then
  git commit todo -m 'Add task'
elif [[ $1 == '-r' ]] | [[ $1 == 'r' ]]; then
  git commit todo -m 'Remove task'
else
  git commit -am "$1"
fi
