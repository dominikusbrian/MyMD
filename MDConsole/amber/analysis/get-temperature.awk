BEGIN{
s=0
}
{
    if ($0 ~ /^ NSTEP/) {
        s = $9
        printf "%f \n", s
        s = 0
    }
}

