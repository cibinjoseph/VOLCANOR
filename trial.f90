program trial
  open(61,file='data.txt',action='write',position='append')
  write(61,*) 'hey'
  close(61)
end program trial
