function device = stimPulse_cpod(~, device)
% send a TTL pulse to a c-pod / cedrus XID-compatible device 
write(device,sprintf("mh%c%c", 255, 0), "char");
end