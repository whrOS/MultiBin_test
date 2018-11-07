function [ec, PSNR, T] = EmbedProc_test(I, a, b, Sele, NL, Tlog, Payload)

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
                if dmax == 0
                    EC(i,j) = EC(i,j) + 1;
                    ED(i,j) = ED(i,j) + 0.5;
                else
                    ED(i,j) = ED(i,j) + 1;
                end
                if dmin == 0
                    EC(i,j) = EC(i,j) + 1;
                    ED(i,j) = ED(i,j) + 0.5;
                else
                    ED(i,j) = ED(i,j) + 1;
                end
            elseif TSele(NL(i,j)+1) == 2
                if dmax == 1
                    EC(i,j) = EC(i,j) + 1;
                    ED(i,j) = ED(i,j) + 0.5;
                else
                    if dmax > 1
                        ED(i,j) = ED(i,j) + 1;
                    end
                end
                if dmin == 1
                    EC(i,j) = EC(i,j) + 1;
                    ED(i,j) = ED(i,j) + 0.5;
                else
                    if dmin > 1
                        ED(i,j) = ED(i,j) + 1;
                    end
                end
            end
        end
    end
end
% sum(EC(:))

flag = 0;
ec = 0;
ed = 0;
for i = 1:floor((A-2)/a)
    if flag == 1
        break
    end
    for j = 1:floor((B-2)/b)
        if NL(i,j) < Tmax
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
T = Tmax;
PSNR = 10*log10(A*B*255^2 / ed);
% ed/ec

end