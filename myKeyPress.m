function myKeyPress(hFig, EventData)
    set(hFig, 'UserData', double(EventData.Key))
end