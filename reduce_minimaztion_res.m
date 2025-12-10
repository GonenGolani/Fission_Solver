function [new_minimaztion_res] = reduce_minimaztion_res(old_minimaztion_res,factor,coeff)

new_minimaztion_res=old_minimaztion_res;
if (coeff==5)
    new_minimaztion_res.MP5=old_minimaztion_res.MP5*factor;
    new_minimaztion_res.proxy5=old_minimaztion_res.proxy5*factor;
    new_minimaztion_res.distal5=old_minimaztion_res.distal5*factor;
end
if (coeff==6)
    new_minimaztion_res.MP6=old_minimaztion_res.MP6*factor;
    new_minimaztion_res.proxy6=old_minimaztion_res.proxy6*factor;
    new_minimaztion_res.distal6=old_minimaztion_res.distal6*factor;
end
if (coeff==7)
    new_minimaztion_res.MP7=old_minimaztion_res.MP7*factor;
    new_minimaztion_res.proxy7=old_minimaztion_res.proxy7*factor;
    new_minimaztion_res.distal7=old_minimaztion_res.distal7*factor;
end
if (coeff==8)
    new_minimaztion_res.MP8=old_minimaztion_res.MP8*factor;
    new_minimaztion_res.proxy8=old_minimaztion_res.proxy8*factor;
    new_minimaztion_res.distal8=old_minimaztion_res.distal8*factor;
end
if (coeff==9)
    new_minimaztion_res.MP8=old_minimaztion_res.MP8*factor;
    new_minimaztion_res.proxy8=old_minimaztion_res.proxy8*factor;
    new_minimaztion_res.distal8=old_minimaztion_res.distal8*factor;
end
if (coeff==10)
    new_minimaztion_res.MP9=old_minimaztion_res.MP9*factor;
    new_minimaztion_res.proxy9=old_minimaztion_res.proxy9*factor;
    new_minimaztion_res.distal9=old_minimaztion_res.distal9*factor;
end
end

