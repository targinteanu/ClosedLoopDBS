function connect_AO()

MAC_addr = 'C8:DF:84:F8:9F:E2'; % MAC address of Functional Neuro Lab's AlphaRS 

Result = AO_DefaultStartConnection(MAC_addr);

if Result
    error(['Connection to AO failed with code ',num2str(Result)])
end

if ~AO_IsConnected()
    error('Did not connect.')
end

end