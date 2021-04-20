for dir in *.case; do
  rm ${dir}/GMCamber
  cd ${dir}; ln -s ../GM21Camber; cd ..
done
