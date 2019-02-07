# MCB432

### Exercise on Feb 5, 2019

Awk commands:

```
awk '/>/{print $1,$2,$3}!/>/{print}' <tmp.txt>tmp2.txt
awk '/>/{print $1,$2,$3}!/>/{print}' <tmp.txt>tmp2.txt
awk '/>/{print $1, substr($2,1,1)".",$3}!/>/{print}' < tmp2.txt > tmp3.txt
```
