function runStim(UserQueue)

StimQueue = parallel.pool.DataQueue;

    QueuedStim = timer(...
        'StartDelay', 10, ...
        'TimerFcn',   {@myPULSE, hObject}, ...
        'StopFcn',    {@finishPULSE, hObject}, ...
        'StartFcn',   {@schedulePULSE, hObject}, ...
        'UserData',   -1);

afterEach(StimQueue, @updateQueuedStim);

% parfeval looped ReadBrain on backgroundPool 

function updateQueuedStim(QueuedData)

% update timer 

% evalin timer to base workspace

end

end