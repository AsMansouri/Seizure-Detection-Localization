function [ Center ] = ChannLoc( Type )
%function [ Center ] = Centers( Type )
%   Give the location of the Electrodes

switch Type
    case 1
        Center = [-145  ,     0;
                   145  ,     0;
                  -85   ,    -270;
                   85   ,    -270;
                 %-413.7333,   -2.9111;
                 %412.8872,   -1.4555;
                   0    ,     0;
                  -115  ,     140;
                   115  ,     140;
                  -230  ,     160;
                   230  ,     160;
                   0    ,     135;
                  -80   ,     260;
                   80   ,     260;
                   0    ,     280;
                  -115  ,    -145;
                   115  ,    -145;
                   0    ,    -135;
                  -280  ,     0;  
                   280  ,     0;
                  -215  ,    -135;
                   215  ,    -135;
                  -285  ,     195;
                   285  ,     195;
                  -253  ,     80;
                   253  ,     80;
                  -350  ,     0;
                   350  ,     0;
                   0    ,    -280];
    case 2
        Center = [  -40,80;
                    -70,50;
                    -80,-25;
                    -50,-65;
                    -25,65;
                    -50,30;
                    -50,-20;
                    -30,-60;
                     35,60;
                     50,25;
                     50,-25;
                     35,-60;
                     50,70;
                     80,40;
                     80,-30;
                     45,-70;
                     0, 40;
                    -0,-40;
                    -100,30;
                     0,100;
                     100,10];
    case 10
        Center = [477.258409785933,120.883027522936;
                  498.146788990826,235.219418960244;
                  522.333333333333,360.549694189602;
                  501.444954128440,484.780581039755;
                  598.191131498471,215.430428134556;
                  645.464831804281,360.549694189602;
                  598.191131498471,507.867737003058;
                  478.357798165138,598.017584097859;
                  314.548929663609,118.684250764526;
                  294.759938837920,235.219418960244;
                  272.772171253823,358.350917431192;
                  293.660550458716,488.078746177370;
                  193.616207951070,218.728593272171;
                  144.143730886850,356.152140672783;
                  195.814984709480,507.867737003058;
                  321.145259938838,600.216360856269];
end
end
