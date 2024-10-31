% Testing: why does myparlfunc keep sending data to the queue (DQ) when it
% is nonempty (i.e. dataQueue.QueueLength = 100 and keeps growing)?

% Start parallel pool if not already started
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool;
end

% Set up a data queue to communicate between the background task and the main thread
dataQueue = parallel.pool.PollableDataQueue;
dataReady = parallel.pool.Constant(false);

% Execute the data acquisition function asynchronously
f = parfeval(pool, @myparlfunc, 0, dataQueue, dataReady);  % 0 indicates no output

iter = 0;
while iter < 100
    pause(.1);
    l = dataQueue.QueueLength
    dataReady.Value = true;
    poll(dataQueue)
    dataReady.Value = false;
    iter = iter+1;
end

cancel(f)

function myparlfunc(DQ, C)
while true
    pause(.01);
    if C.Value
        send(DQ, rand);
    end
end
end