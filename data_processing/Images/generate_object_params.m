clear all
close all

cd '/Users/James/James_scripts/data_processing/Images'

%[white_bufer gamma_correction boundary_thickness object_topology]
obj_params = ...
    [10 1 2 0; % apple.jpg
    10 0.75 1 1;% armadillo.jpg
    10 0.75 1 1;% artichoke.jpg
    10 0.75 1 1;% asparagus.jpg
    10 0.75 1 1;% banana.jpg
    10 0.75 1 1;% barn_owl.jpg
    10 0.75 1 1;% bat.jpg
    10 0.75 1 1;% bee.jpg
    10 0.6 1 0;% beetle.jpg
    10 0.75 1 1;% bellflowers.jpg
    10 0.75 1 1;% butterfly.jpg
    10 0.75 1 0;% cabbage.jpg
    10 0.75 1 1;% cat.jpg
    10 0.75 2 1;% cedar.jpg
    10 0.625 2 1;% celery.jpg
    10 0.75 1 0;% chard.jpg
    10 1 1 1;% cheetah.jpg
    10 0.75 1 0;% cockle.jpg
    10 0.5 1 1;% cow.jpg
    10 0.75 1 1;% crocodile.jpg
    10 0.75 1 0;% cucumber.jpg
    10 0.75 1 1;% dolphin.jpg
    10 0.55 2 0;% dromedary.jpg
    10 0.5 1 1;% duck.jpg
    5 0.5 2 0;% eggplant.jpg
    10 0.75 1 0;% elephant.jpg
    10 0.7 3 1;% giraffe.jpg
    10 0.75 1 1;% goldfinch.jpg
    0 0.5 3 0;% goose.jpg
    1 1 2 1;% grapes.jpg
    5 0.8 1 1;% hen.jpg
    10 0.5 1 1;% hippopotamus.jpg
    3 0.75 1 0;% horse.jpg
    5 0.5 1 0;% hummingbird.jpg
    10 0.5 1 0;% kangaroo.jpg
    10 0.75 1 1;% kille_ whale.jpg
    10 0.75 2 1;% kiwi.jpg
    10 0.5 1 0;% leek.jpg
    10 0.7 2 0;% lemon.jpg
    10 0.75 1 0;% lettuce.jpg
    3 0.75 1 0;% lioness.jpg
%     1 0.5 2 1;% lobster.jpg
    10 0.5 1 1;% lynx.jpg
    10 0.4 1 1;% magpie.jpg
    10 0.75 1 1;% manatee.jpg
    10 0.75 1 0;% melon.jpg
    10 0.75 1 0;% narwhal.jpg
    3 0.5 1 0;% onion.jpg
    10 0.75 1 0;% orange.jpg
    10 0.5 1 1;% ostrich.jpg
    10 0.6 1 1;% owl.jpg
    10 0.75 1 0;% pansy.jpg
    10 0.75 1 1;% partridge.jpg
    10 0.75 1 1;% pelican.jpg
    10 0.5 1 0;% penguin.jpg
    5 0.75 2 0;% pepper.jpg
    10 0.5 1 0;% pheasant.jpg
    10 0.5 1 1;% pigeon.jpg
    10 0.5 1 1;% platypus.jpg
    3 0.75 2 0;% pomegranate.jpg
    10 0.5 2 1;% pomfret.jpg
    10 0.5 1 1;% pumpkin.jpg
    10 0.5 2 1;% quince.jpg
    4 0.5 1 1;% redcurrant.jpg
    2 0.5 5 0;% rhino.jpg
    2 0.5 1 1;% rooster.jpg
    10 0.75 1 1;% rose.jpg
    10 0.75 1 1;% seagull.jpg
    10 0.5 1 1;% shark.jpg
    5 0.6 2 1;% snake.jpg
    10 0.5 1 1;% sparrow.jpg
    10 0.5 1 0;% strawberry.jpg
    10 0.75 1 1;% sunflower.jpg
    10 0.5 1 1;% tapir.jpg
    1 0.75 1 1;% tiger.jpg
    10 0.5 1 1;% toucan.jpg
    10 0.6 1 1;% turtle.jpg
    3 0.75 3 0;% watermelon.jpg
    10 0.5 1 1;% whale.jpg
    10 0.6 1 1];% zebra.jpg

save obj_params