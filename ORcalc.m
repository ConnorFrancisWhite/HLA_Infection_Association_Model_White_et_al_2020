function OR = ORcalc(p1,p2);

OR = (p1/(1-p1))/(p2/(1-p2));
if p1 == 0
    OR = 0; 
end
if p2 == 0
    OR = 0; 
end
