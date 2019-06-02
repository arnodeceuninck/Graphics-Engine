for i in *.ini; do
    ../engine $i || break
done
