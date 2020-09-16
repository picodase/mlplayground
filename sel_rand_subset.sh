mkdir subset

a=( * )
randf=( "${a[RANDOM%${#a[@]}]"{1..3}"}" )

for (( i = 0; i <= ${#randf[@]}-1; i++ ))
do
    cp "${randf[i]}" subset
    #printf "${randf[i]}"
done