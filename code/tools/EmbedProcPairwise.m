function [ec, PSNR, T] = EmbedProcPairwise(I, a, b, Sele, NL, Tlog, Payload)
optMapEC0 = zeros(256,256); optMapEC0(1,:) = 1;  optMapEC0(:,1) = 1;  optMapEC0(1,1) = log2(3); optMapEC0(2,2) = 1;
optMapED0 = ones(256,256)*2;optMapED0(1,:) = 3/2;optMapED0(:,1) = 3/2;optMapED0(1,1) = 2/3;     optMapED0(2,2) = 1;

optMapEC1 = zeros(256,256); optMapEC1(2:end, 2:end) = optMapEC0(1:255,1:255); optMapEC1(2,1) = 1;  optMapEC1(1,2) = 1;
optMapED1 = zeros(256,256); optMapED1(2:end, 2:end) = optMapED0(1:255,1:255);
optMapED1(3:end,1) = 1; optMapED1(1,3:end) = 1; optMapED1(2,1) = 0.5; optMapED1(1,2) = 0.5;

[A, B] = size(I);

Tmax = Tlog(find(Sele~=0, 1, 'last' ));

TSele = zeros(1, Tmax);
Cnt = 1;
for i = 1 : 1 : Tmax
    TSele(i) = Sele(Cnt);
    if i >= Tlog(Cnt)
        Cnt = Cnt + 1;
    end
end

EC = zeros(floor((A-2)/a),floor((B-2)/b));
ED = zeros(floor((A-2)/a),floor((B-2)/b));
for i = 1:floor((A-2)/a)
    for j = 1:floor((B-2)/b)
        if NL(i,j) < Tmax
            X = I(a*(i-1)+1:a*i,b*(j-1)+1:b*j);
            X = X(:);
            [Y, In] = sort(X);
            % max
            if In(a*b) < In(a*b-1)
                dmax = Y(a*b) - Y(a*b-1) - 1;
            else
                dmax = Y(a*b) - Y(a*b-1);
            end
            % min
            if In(2) < In(1)
                dmin = Y(2) - Y(1) - 1;
            else
                dmin = Y(2) - Y(1);
            end
            if TSele(NL(i,j)+1) == 1
                EC(i,j) = EC(i,j) + optMapEC0(dmax+1,dmin+1);
                ED(i,j) = ED(i,j) + optMapED0(dmax+1,dmin+1);
            elseif TSele(NL(i,j)+1) == 2
                EC(i,j) = EC(i,j) + optMapEC1(dmax+1,dmin+1);
                ED(i,j) = ED(i,j) + optMapED1(dmax+1,dmin+1);
            end
        end
    end
end


flag = 0;
for T = 1 : 1 : Tmax
    ec = 0;
    ed = 0;
    for i = 1:floor((A-2)/a)
        if flag == 1
            break
        end
        for j = 1:floor((B-2)/b)
            if NL(i,j) < T
                ec = ec + EC(i,j);
                ed = ed + ED(i,j);
                if ec >= Payload
                    Kend = [i,j];
                    flag = 1;
                    break;
                end
            end
        end
    end
    if flag == 1
        break
    end
end

PSNR = 10*log10(A*B*255^2 / ed);
% ed/ec

end