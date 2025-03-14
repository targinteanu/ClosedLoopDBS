function disconnect_AO()

Results = AO_CloseConnection();
if (Results == 0)
    disp('Connection closed successfully.');
else
    % throw error?
    disp(['Connection close error code ',num2str(Results)]);
end

if AO_isConnected()
    error('Did not disconnect.')
end

end