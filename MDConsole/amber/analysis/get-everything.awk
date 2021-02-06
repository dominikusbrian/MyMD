BEGIN{
s=0
}
{
    if ($0 ~ /^ NSTEP/) {
        s = $9
    }
    if ($0 ~ /^ Etot/){
	a = $3
	b = $6 
	c = $9
	printf "%f,%f,%f,%f\n", s, a, b, c
	s = 0 
	a = 0
	b = 0
	c = 0
	}
}

