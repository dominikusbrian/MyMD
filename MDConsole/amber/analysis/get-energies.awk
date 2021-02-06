BEGIN{
s=0
}
{
    if ($0 ~ /^ BOND/) {
        s += $3 + $6 + $9
    }
    if ($0 ~ /^ VDWAALS/) {
        s += $3 + $6 + $9
    }
    if ($0 ~ /^ 1-4 VDW/) {
        s += $4 + $8 + $11
        printf "%f \n", s
        s = 0
    }
}
