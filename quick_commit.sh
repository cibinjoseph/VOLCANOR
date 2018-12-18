case $1 in
  -a|--add)
    git commit todo -m "Add task" ;;
  -h|--help)
    echo "Usage: quick_commit {--add|--remove|--test}" ;;
  -r|--remove)
    git commit todo -m "Remove task" ;;
  -t|--test)
    git commit -am "testcase" ;;
  *)
    if [ $# -eq 0 ]; then  # No arguments
      git commit -am "Quick commit"
    else
      git commit -am "$1"
    fi
esac
