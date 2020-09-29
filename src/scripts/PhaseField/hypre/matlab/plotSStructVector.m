function plotSStructVector(x,rfactors,labels)

ndim = size(rfactors,2);

switch ndim
    case 2
        i=1;
        rank0=1;
        r = [1 1];

        for l=1:length(x.StructVector)
            r = r .* rfactors(l,1:2);
            for p=1:x.StructVector(l).NParts
                sol(i).Low = x.StructVector(l).Parts(p).Low;
                sol(i).High = x.StructVector(l).Parts(p).High;
                sol(i).Level = l;
                dim = sol(i).High-sol(i).Low+3;
                rank1 = rank0+prod(dim)-1;
                sol(i).X = ([sol(i).Low(1):sol(i).High(1)]+0.5)/r(1);
                sol(i).Y = ([sol(i).Low(2):sol(i).High(2)]+0.5)/r(2);
                c = reshape(x.Full(rank0:rank1),dim')';
                c([1 end],:) = [];
                c(:,[1 end]) = [];
                sol(i).C = c;
                lbl = reshape(labels(rank0:rank1),dim');
                lbl([1 end],:) = [];
                lbl(:,[1 end]) = [];
                sol(i).Labels = string(reshape(lbl,[],1));
                i=i+1;
                rank0=rank1+1;
            end
        end

        for l=1:length(x.StructVector)
            hold on
            m=find([sol.Level]==l);
            for i=1:length(m)
                image(sol(m(i)).X,sol(m(i)).Y,sol(m(i)).C,'CDataMapping','scaled')
                [X,Y]=ndgrid(sol(m(i)).X,sol(m(i)).Y);
                text(reshape(X,[],1),reshape(Y,[],1),sol(m(i)).Labels,'HorizontalAlignment','center','VerticalAlignment','middle');
            end
        end

        xlabel('x')
        ylabel('y')

    case 3
        i=1;
        rank0=1;
        r = [1 1 1];

        T = [0;0;0];
        for p=1:x.StructVector(1).NParts
            T = max(T,(x.StructVector(1).Parts(p).High + 1) ./ (2*r(1,:))');
        end

        ax = gca;
        cols = ax.ColorOrder;
        ncols = size(cols,2);

        for l=1:length(x.StructVector)
            r = r .* rfactors(l,:);
            color = cols(mod(l,ncols)+1,:);

            for p=1:x.StructVector(l).NParts
                sol(i).Low = x.StructVector(l).Parts(p).Low;
                sol(i).High = x.StructVector(l).Parts(p).High;
                sol(i).Level = l;
                dim = sol(i).High-sol(i).Low+3;
                rank1 = rank0+prod(dim)-1;

                sol(i).X = ([sol(i).Low(1):sol(i).High(1)+1])/r(1) - T(1);
                sol(i).Y = ([sol(i).Low(2):sol(i).High(2)+1])/r(2) - T(2);
                sol(i).Z = ([sol(i).Low(3):sol(i).High(3)+1])/r(3) - T(3);

                plot3( ...
                    sol(i).X([ 1   1   1   1   1 ]), ...
                    sol(i).Y([ 1   1  end end  1 ]), ...
                    sol(i).Z([ 1  end end  1   1 ]), ...
                    'Color', color )
                hold on

                plot3( ...
                    sol(i).X([end end end end end]), ...
                    sol(i).Y([ 1   1  end end  1 ]), ...
                    sol(i).Z([ 1  end end  1   1 ]), ...
                    'Color', color )

                plot3( ...
                    sol(i).X([ 1  end]), ...
                    sol(i).Y([ 1   1 ]), ...
                    sol(i).Z([ 1   1 ]), ...
                    'Color', color )

                plot3( ...
                    sol(i).X([ 1  end]), ...
                    sol(i).Y([ 1   1 ]), ...
                    sol(i).Z([end end]), ...
                    'Color', color )

                plot3( ...
                    sol(i).X([ 1  end]), ...
                    sol(i).Y([end end]), ...
                    sol(i).Z([ 1   1 ]), ...
                    'Color', color )

                plot3( ...
                    sol(i).X([ 1  end]), ...
                    sol(i).Y([end end]), ...
                    sol(i).Z([end end]), ...
                    'Color', color )

                c = reshape(x.Full(rank0:rank1),dim');
                c([1 end],:,:) = [];
                c(:,[1 end],:) = [];
                c(:,:,[1 end]) = [];
                sol(i).C = c;
                lbl = reshape(labels(rank0:rank1),dim');
                lbl([1 end],:,:) = [];
                lbl(:,[1 end],:) = [];
                lbl(:,:,[1 end]) = [];
                sol(i).Labels = lbl;
                i=i+1;
                rank0=rank1+1;
            end
        end

        for l=1:length(x.StructVector)
            m=find([sol.Level]==l);
            for i=1:length(m)
                x0 = find( sol(m(i)).X >= 0, 1);
                y0 = find( sol(m(i)).Y >= 0, 1);
                z0 = find( sol(m(i)).Z >= 0, 1);
                X = ( sol(m(i)).X(1:end-1)+sol(m(i)).X(2:end) )/2;
                Y = ( sol(m(i)).Y(1:end-1)+sol(m(i)).Y(2:end) )/2;
                Z = ( sol(m(i)).Z(1:end-1)+sol(m(i)).Z(2:end) )/2;
                [X,Y,Z]=ndgrid(X,Y,Z);
                if x0
                    xImage = [0,0;0,0];
                    yImage = [sol(m(i)).Y([1 1]);sol(m(i)).Y([end end])];
                    zImage = [sol(m(i)).Z([1 end]);sol(m(i)).Z([1 end])];
                    cImage = squeeze(sol(m(i)).C(x0,:,:));
                    surf(xImage,yImage,zImage,...    %# Plot the surface
                        'CData',cImage,...
                        'EdgeColor',color,...
                        'FaceAlpha',0.5,...
                        'FaceColor','texturemap');

                    text(...
                        reshape(0*X(x0,:,:),[],1),...
                        reshape(Y(x0,:,:),[],1),...
                        reshape(Z(x0,:,:),[],1),...
                        string(reshape(sol(m(i)).Labels(x0,:,:),[],1)),...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','middle'...
                        );
                end
                if y0
                    xImage = [sol(m(i)).X([1 1]);sol(m(i)).X([end end])];
                    yImage = [0,0;0,0];
                    zImage = [sol(m(i)).Z([1 end]);sol(m(i)).Z([1 end])];
                    cImage = squeeze(sol(m(i)).C(:,y0,:));
                    surf(xImage,yImage,zImage,...    %# Plot the surface
                        'CData',cImage,...
                        'EdgeColor',color,...
                        'FaceAlpha',0.5,...
                        'FaceColor','texturemap');

                    text(...
                        reshape(X(:,y0,:),[],1),...
                        reshape(0*Y(:,y0,:),[],1),...
                        reshape(Z(:,y0,:),[],1),...
                        string(reshape(sol(m(i)).Labels(:,y0,:),[],1)),...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','middle'...
                        );
                end
                if z0
                    xImage = [sol(m(i)).X([1 1]);sol(m(i)).X([end end])];
                    yImage = [sol(m(i)).Y([1 end]);sol(m(i)).Y([1 end])];
                    zImage = [0,0;0,0];
                    cImage = squeeze(sol(m(i)).C(:,:,z0));
                    surf(xImage,yImage,zImage,...    %# Plot the surface
                        'CData',cImage,...
                        'EdgeColor',color,...
                        'FaceAlpha',0.5,...
                        'FaceColor','texturemap');

                    text(...
                        reshape(X(:,:,z0),[],1),...
                        reshape(Y(:,:,z0),[],1),...
                        reshape(0*Z(:,:,z0),[],1),...
                        string(reshape(sol(m(i)).Labels(:,:,z0),[],1)),...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','middle'...
                        );
                end
            end
        end

        xlabel('x')
        ylabel('y')
        zlabel('z')

end
end