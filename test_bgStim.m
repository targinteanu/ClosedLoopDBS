function T = test_bgStim()
stimulator = defineSTIM4(1,2,1000,1000,100,100,100,50,1);
tic;
stimulator.play(1);
T = toc;
stimulator.disconnect;
end