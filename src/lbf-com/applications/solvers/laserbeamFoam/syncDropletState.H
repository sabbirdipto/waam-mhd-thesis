// Persist the droplet accumulator so a restart does not discard it.
// Included immediately after updateDropletSource.H, so the values written are
// the ones that step just produced. See patch_droplet_restart.py.
{
    if (dropOn)
    {
        dropletState.set("dropOwed",      dropOwed_);
        dropletState.set("dropImpactEnd", dropImpactEnd_);
        dropletState.set("dropNum",       dropNum_);
    }
}
